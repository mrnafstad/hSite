<template>
	<div>
		<div id="blog">
			Blogg
			<div v-if="allPosts">
				<li v-for="(post, idx) in allPosts" :key="idx">
					<h3>{{ post.title }}</h3>
					<p>{{ post.timeOfPost.toDate().toDateString() }}</p>
					<p>{{ post.post }}</p>
					<!--<button id="likes" @click="like(post)">Like! ({{ post.likes }})</button>-->
				</li>
			</div>
		</div>
		<NewPost v-if="auth" />
	</div>
</template>

<script>
import { mapGetters,mapActions } from 'vuex'
import NewPost from '../components/NewPost.vue'

export default {
	name: 'Blog',
	components: {
		NewPost
	},
	methods: {
		...mapActions(['getPosts', 'newLike']),
		like(post) {
			
			console.log(post)
			this.newLike(post.id)
		}
	},
	computed: {
		...mapGetters(['allPosts', 'auth'])
	},
	created() {
		this.getPosts().then(console.log("posts fetched"))
	}
}
</script>

<style>
	#blog{
		border-style: outset;
		background-color: snow;
	}
	#likes {
		border-style: outset;
		background-color: aqua;
		border-color: lightblue;
	}
</style>